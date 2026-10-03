# Digest — Zhang, Li, Baoyin (2025), "Noncoplanar Capture Trajectory Optimization in Uranian System Exploration Mission" (IEEE TAES)

**Digested:** 2026-10-03 (text-layer PDF, 11 pages, read in full; the two-column text layer was cross-read against
the rendered pages; Tables I-IV, Figures 6-9 were read from the rendered pages because they are images).
**Purpose:** record exactly what in-system trajectory content at Uranus is published, ahead of work on the
catalogued Uranian two-moon quasi-cyclers. Facts only; no novelty verdict is written here.

## Citation
Nan Zhang, Haiyang Li, Hexi Baoyin, "Noncoplanar Capture Trajectory Optimization in Uranian System Exploration
Mission", IEEE Transactions on Aerospace and Electronic Systems, vol. 61, no. 6, December 2025, pp. 18732-18742
(running header: "IEEE TRANSACTIONS ON AEROSPACE AND ELECTRONIC SYSTEMS VOL. 61, NO. 6 DECEMBER 2025").

- Author names printed in capitals: NAN ZHANG (Tsinghua University, Beijing, China); HAIYANG LI (Deep Space
  Exploration Laboratory, "Heifei", China, as printed); HEXI BAOYIN (Tsinghua University, Beijing, China).
  Corresponding author: Baoyin.
- Dates as printed: "Received 29 April 2025; revised 12 August 2025; accepted 22 September 2025. Date of
  publication 7 October 2025; date of current version 8 December 2025."
- DOI as printed: "DOI. No. 10.1109/TAES.2025.3618517". ISSN line: "0018-9251 (c) 2025 IEEE". "Refereeing of this
  contribution was handled by M. A. Ayoubi."
- Funding as printed: National Natural Science Foundation of China Grant 12402052; Postdoctoral Fellowship
  Program of CPSF Grant BX20240201; Frontier Research Program of the Deep Space Exploration Laboratory Grant
  GC04FY1003ZC3ZT-2326. Acknowledgement: the authors acknowledge ChatGPT for language editing.
- The PDF pages carry an IEEE Xplore download footer (institutional licence stamp, downloaded 4 August 2026);
  this is a footer, not part of the paper.
- Filed in the private paper corpus as
  `zhang-li-baoyin-2025-noncoplanar-capture-trajectory-optimization-uranian-system-exploration-ieee-taes-doi-10.1109-TAES.2025.3618517.pdf`,
  md5 `80f5912bdd0d30903a6cf6f7e0180d7e`.

## What the paper is
Abstract (quoted): "The noncoplanar capture trajectory optimization problem involves guiding a spacecraft into a
designated noncoplanar orbital plane, such as the equatorial plane, starting from a given incoming excess
velocity vector. This challenge is particularly significant for the Uranian system exploration mission, due to
the near-perpendicular orientation of Uranian equatorial plane relative to the ecliptic plane. This article
formulates the problem and investigates two capture strategies: the multiimpulse transfer strategy and the
Titania resonant gravity-assist (GA) strategy. Global and local continuous variable optimization methods are
employed to solve the impulse trajectory optimization problem, while the A* algorithm is used to address the
resonant GA trajectory optimization problem. Simulation results indicate that the three-impulse transfer
strategy is more straightforward to implement and does not incur excessive Delta-v costs, making it a more
suitable option when the initial angle is small and the transfer time is short. Conversely, the Titania
resonant GA strategy offers greater Delta-v savings, making it potentially the better choice when the initial
angle is large and the transfer time is long."

Sections as printed: I Introduction; II Origin of Noncoplanar Capture Issue; III Noncoplanar Capture Trajectory
Optimization Problem; IV Two Noncoplanar Capture Strategies (A Multiple Impulse Transfer Strategy; B Moon
Resonant GA Transfer Strategy); V Simulation Results (A Multiple Impulse Transfer Trajectory; B Titania
Resonant GA Trajectory; C Discuss on the Influence of Excess Velocity Vector); VI Conclusion. Four tables (I-IV),
nine figures, 33 references.

## What "capture" means here
Capture is the transition "from an incoming hyperbolic trajectory to a bound elliptical orbit around the target
planet" (Introduction), considered in this paper as a Uranus-centred problem that starts at the Uranian sphere of
influence and ends in a bound 15-day equatorial orbit. Quoted definition: the spacecraft "at the initial time
t0 ... is on a hyperbolic trajectory with an incoming excess vector v-inf, just arriving at the boundary of the
planetary sphere of influence at a distance r0 from the planet's center. Through a series of maneuvering
operations, the spacecraft is eventually required to enter a target elliptical orbit that lies in a noncoplanar
plane from the initial orbital plane." Target orbit parameters (Table I): a_f = 6.2701e5 km, e_f = 0.9553;
text: "The parameters a_f and e_f of the target orbit are chosen such that its periapsis exactly matches r_min
and the orbital period is 15 days." The target plane is Uranus's equatorial plane (i_f = 0, RAAN_f = 0 in the
chosen frame, Eq. 10).

Only Uranus orbit insertion by impulsive maneuvers is optimised, plus, in one strategy, Titania gravity assists
interleaved with the impulses. "Moon-aided capture" is Titania only: "The GA effect of other moons is weaker than
that of Titania, as the second-largest moon, Oberon, although similar in mass, has a longer orbital period,
requiring more time for resonant GAs, while the remaining moons have masses far smaller than Titania's."
Aerocapture, low-thrust capture and ballistic capture (weak stability boundaries, invariant manifolds, cited for
Mars) are discussed in the introduction and not used.

The UOP context, quoted: "To date, the only proposed noncoplanar capture case appears in the UOP mission
concept. It plans an initial impulsive capture maneuver, followed by a series of gravity assists (GAs) from the
Uranian moon Titania to achieve plane alignment over nearly two years [3]."

## Model
- Two-body (Uranus) Keplerian arcs: "ṙ = v, v̇ = -(mu/r^3) r" (Eq. 1); "The orbit of the unpowered spacecraft is
  only affected by the gravitational force of the central planet". Impulsive maneuvers (Eq. 2).
- Moon GAs: "The linked-conics model is adopted to simplify the computation of GAs: the spacecraft's velocity
  vector is assumed to change instantaneously, while its position and the magnitude of its hyperbolic excess
  velocity relative to the moon remain unchanged" (Eq. 14); bending limited by delta_max = 2 arcsin(mu_m /
  (mu_m + r_m,min v_inf,m^2)) (Eq. 15). Titania is treated as on a circular equatorial orbit.
- Heliocentric stage (section II only): planetary positions from DE430; two-body approximation; linked-conics
  (zero-sphere-of-influence patched-conics) GAs; C3 matching method [28] to generate trajectories exhaustively.
- Optimisation: the multiple-impulse problem uses Keplerian propagation plus a Lambert arc; particle swarm
  optimisation (PSO) for the initial guess then the SBPLX local algorithm (NLopt); for the resonant GA problem an
  "A* algorithm with pruning [24]" is the subroutine for the multiple resonant GAs, nested inside PSO and SBPLX.
  No CR3BP, no Tisserand graph, no invariant manifolds, no pump-and-crank v-infinity-sphere method.

## Section II: arrival geometry at Uranus
- The moons orbit in the equatorial plane: "the Jovian Galilean moons, the Saturnian Titan and Enceladus, and
  the Uranian Ariel, Umbriel, Titania, and Oberon all follow this pattern."
- Four heliocentric GA sequences are scanned: EVEEJ (Jupiter, Galileo-like), EVVEJS (Saturn, Cassini-like),
  EVEEJU and E-v-EJU (Uranus; "reference the UOP mission concept [3]"); launch 2028-2032, launch C3 10-20 km2/s2
  (E-v-EJU 20-30), 30-day launch discretisation, 0.2 km2/s2 C3 step; duration caps 7 years (Jupiter), 12 years
  (Saturn), 16 years (Uranus); delta-v cap 1000 m/s for E-v-EJU.
- Result: for Uranus missions the angle theta between v-inf and the equatorial plane "can reach up to 50 deg"
  (Fig. 1, all recorded trajectories) and up to "as high as 80 deg" over 100 years of arrival times (Fig. 2); the
  Uranian equatorial plane is "inclined at an extreme angle of 82.7 deg to the ecliptic"; theta varies with a
  period of 84.3 years (trigonometric fit, amplitude 82.7 deg).
- Example heliocentric trajectory (Fig. 3): "Earth-Venus-Earth-Earth-Jupiter-Uranus", departed 19 November 2029,
  15.7-year journey, arrival v-inf of about 6 km/s and theta of about 60 deg; this is the source of the
  v-inf = 6 km/s, theta = 60 deg test case.

## Section V: results

### Table I, exact transcription ("Uranus Related Parameters")
| parameter | mu | r0 | a_f | e_f | r_min | mu_m | a_m | r_m,min |
|---|---|---|---|---|---|---|---|---|
| value | 5.7940e6 | 5.1743e7 | 6.2701e5 | 0.9553 | 2.800e4 | 226.90 | 4.3630e5 | 838.9 |
| unit | km^3/s^2 | km | km | / | km | km^3/s^2 | km | km |

Text: "mu_m, a_m, and r_m,min denote the gravitational constant, the semimajor axis, and the minimum flyby radius
of Titania, respectively" (the last three table columns are the Titania entries; mu, r0, a_f, e_f and r_min are
Uranus and target-orbit entries).

### Test case for Tables II-IV: v-inf = 6 km/s, theta = 60 deg, t_max = 360 days
Coplanar reference, Eq. (17): "Delta-v = sqrt(2mu/r_min + v_inf^2) - sqrt(mu(1 + e_f)/r_min) = 1093.8 m/s".

**Table II, "Details of Three-Impulse Noncoplanar Capture Trajectory"**
| Event | Time (days) | Inclination after impulse | Delta-v (m/s) |
|---|---|---|---|
| Initial value | 0 | 86.5 deg | 0 |
| First impulse | 98.2 | 86.5 deg | 1128.4 |
| Second impulse | 122.1 | 0 deg | 206.3 |
| Third impulse | 360.0 | 0 deg | 194.6 |
| Total | - | - | 1529.3 |

**Table III, "Details of Four-Impulse Noncoplanar Capture Trajectory"**
| Event | Time (days) | Inclination after impulse | Delta-v (m/s) |
|---|---|---|---|
| Initial value | 0 | 87.6 deg | 0 |
| First impulse | 98.2 | 87.6 deg | 922.6 |
| Second impulse | 126.0 | 49.9 deg | 240.5 |
| Third impulse | 330.6 | 0 deg | 106.1 |
| Four impulse | 360.0 | 0 deg | 194.4 |
| Total | - | - | 1463.6 |

**Table IV, "Details of Titania GA Noncoplanar Capture Trajectory"**
| Event | Time (days) | Inclination after impulse/GA | Delta-v (m/s) |
|---|---|---|---|
| Initial value | 0 | 69.9 deg | 0 |
| First impulse | 98.2 | 69.9 deg | 1034.1 |
| Second impulse | 109.1 | 35.7 deg | 245.4 |
| Titania GA (3:1) | 126.6 | 33.4 deg | 0 |
| Titania GA (3:1) | 152.7 | 30.0 deg | 0 |
| Titania GA (3:1) | 178.8 | 26.4 deg | 0 |
| Titania GA (3:1) | 204.9 | 22.5 deg | 0 |
| Titania GA (3:1) | 231.1 | 18.3 deg | 0 |
| Titania GA (3:1) | 257.2 | 13.9 deg | 0 |
| Titania GA (3:1) | 283.3 | 9.3 deg | 0 |
| Titania GA (3:1) | 309.4 | 4.7 deg | 0 |
| Titania GA | 335.5 | 0 deg | 0 |
| Third impulse | 347.4 | 0 deg | 0.005 |
| Four impulse | 360.0 | 0 deg | 67.2 |
| Total | - | - | 1346.8 |

Arithmetic note (reader's check, not in the paper): the Table IV Delta-v entries sum to 1346.705 m/s against the
printed total 1346.8; the Table II and III entries sum exactly to their printed totals. Consecutive Titania GA
times in Table IV are spaced 26.1 days (26.2 once).

Text for the Titania case, quoted: "after the initial capture maneuver, the spacecraft performs a second impulse
to change the orbital inclination while transferring toward Titania, arriving 126.6 days later. Subsequently,
nine Titania GAs are used to gradually align the orbital plane with Uranian equatorial plane. The first eight
GAs place the spacecraft into a 3:1 resonance orbit with Titania. During each Titania flyby, the relative
velocity is about 4.885 km/s, and the maximum velocity change achievable through GA is 109.5 m/s. Finally, two
additional impulsive maneuvers are executed to insert the spacecraft into the target orbit. Compared with the
optimal multiple-impulse transfer trajectory, the Titania GA transfer trajectory further reduces the total
Delta-v from 1463.6 to 1346.8 m/s. This reduction is primarily attributed to the multiple Titania flybys, which
significantly decreased the orbital inclination." The paper does not define which side of "3:1" is the
spacecraft period; Figure 8(d) plots the period over the resonant-GA interval and Figure 8(b) shows the
distance to Uranus (about 2e6 km maximum) with nine GA markers at Titania's distance.

Other quoted statements for the three-impulse case: the first impulse "is performed 34 deg past periapsis" and
contributes most of the Delta-v; "the optimal initial inclination is not equal to theta (= 60 deg), but rather a
significantly larger angle of 86.5 deg"; the noncoplanar penalty versus coplanar is "an additional Delta-v cost
of 435.5 m/s, along with an additional transfer time of 261.8 days (calculated as 360 - 98.2 days)". Four-
impulse first impulse "23 deg past periapsis"; Titania case first impulse "17 deg past periapsis". The paper
also reports that two-impulse and five-impulse trajectories found no solution outperforming the four-impulse one.

### Parametric study (Figure 9, nine panels; no table)
v-inf = 4, 6, 8 km/s; theta = 10, 50, 90 deg; t_max 150 to 1000 days (step 10 days for the multiple-impulse
strategies, 100 days for the Titania GA strategy). Panel layout, as read from the text references:
(a), (d), (g) are theta = 10 deg; (b), (e), (h) theta = 50 deg; (c), (f), (i) theta = 90 deg; rows are v-inf = 4,
6, 8 km/s. Statements: at theta = 10 deg the three strategies are not significantly different from coplanar
capture, and "When tf > 500 days, the Titania GA transfer may yield a slightly lower Delta-v than the coplanar
capture, but the advantage remains marginal"; for theta = 50 and 90 deg Delta-v increases significantly and
decreases with longer durations, "with the Titania GA transfer offering notable Delta-v savings"; increasing
v-inf "primarily raises the required Delta-v". Values read by eye from Figure 9 (approximate, not tabulated in
the paper): for v-inf = 8 km/s, theta = 90 deg (panel i) the three-impulse curve is near 3100 m/s while the
Titania GA points fall from about 3200 m/s at 200 days to about 1900 m/s at 1000 days, with the coplanar line
near 1750 m/s.

Conclusion, quoted: "Simulation results show that when the initial angle is small and the transfer time is
short, the required Delta-v differences among three strategies are small, making the three-impulse strategy more
favorable due to its engineering simplicity. Conversely, under conditions of larger angle and longer transfer
duration, the resonant GA strategy demonstrates greater Delta-v savings and may be the more optimal choice."
Computation cost: a multiple-impulse case under 0.5 s; a resonant-GA case up to 2 minutes (i7-7700), large
sweeps on a supercomputing platform.

## V-infinity at moons and resonances (summary)
- Titania only; relative velocity at each Titania flyby about 4.885 km/s; maximum GA velocity change 109.5 m/s
  (a per-flyby bound; reader's check, not in the paper: Eq. 15 with the Table I Titania values and
  v-inf,m = 4.885 km/s gives delta_max of about 0.0224 rad and 2 v-inf sin(delta_max/2) of about 109.5 m/s,
  matching the printed figure).
- Resonance: "3:1" (eight GAs labelled, ninth unlabelled in Table IV); no other resonance appears; Oberon,
  Umbriel, Ariel are not flown by.
- Durations: first Titania encounter at 126.6 days, last Titania GA at 335.5 days, final orbit insertion at 360.0
  days (t_max); Titania sequence spans 209 days (335.5 - 126.6, reader's arithmetic).
- Delta-v totals for v-inf = 6 km/s, theta = 60 deg, t_max = 360 d: coplanar 1093.8; three-impulse 1529.3;
  four-impulse 1463.6; Titania GA 1346.8 m/s.

## Term search (grep -a over layout text, line-joined, ligatures normalised)
cycler 0; cycling 0; periodic 0; repeating 0; free-return 0; resonant 22; resonance 2; heteroclinic 0;
homoclinic 0; quasi-periodic 0; torus 0; four-body 0; Umbriel 1; Titania 30; Oberon 2; Ariel 1; Miranda 0;
Triton 0; Proteus 0; Nereid 0. (Counts are from the two-column layout text, and include captions, references
and headers.) The single Umbriel and Ariel hits are the list "the Uranian Ariel, Umbriel, Titania, and Oberon all
follow this pattern"; Oberon appears in that list and in the second-largest-moon sentence quoted above.

Sentences for the requested terms: cycler, cycling, repeating, free-return: none, those words do not occur.
The word "resonant" refers to the "Titania resonant GA strategy", i.e. a sequence of resonant gravity assists of
one moon by one spacecraft ("this section explores the use of these moons for multiple resonant GAs to achieve an
inclination change during the capture phase"); "multiple resonant GAs using a single Jovian moon for a 90 deg
inclination change" is attributed to Huang et al. [25].

## Does the paper contain a trajectory that returns to the same two moons, or a periodic or quasi-periodic orbit that encounters two moons
No. One moon (Titania) is flown by nine times in a single 360-day capture sequence that ends in a bound 15-day
orbit at Uranus; no second moon is encountered; no periodic, quasi-periodic or repeating trajectory is
constructed.

## What it does NOT contain
- Any flyby of Ariel, Umbriel, Oberon or Miranda (Oberon is dismissed in one sentence as needing more time).
- Any post-capture moon-tour design beyond the nine-flyby Titania sequence; no tour of the major moons; no
  multi-moon V-infinity matching.
- Any periodic, quasi-periodic, halo, Lyapunov, torus, homoclinic or heteroclinic orbit, any CR3BP or invariant
  manifold model, any four-body model.
- Any V-infinity table at the moons other than the single value "about 4.885 km/s".
- Any Tisserand graph or pump-and-crank diagram; resonance is "3:1" only.
- A UOP-style published trajectory list; the UOP concept is cited as [3] for the "Titania GA over nearly two
  years" idea, and this paper's own capture is 360 days (the paper does not compare the two durations).
- Ephemeris-level (DE430) moon motion for the capture; Titania is circular and equatorial.

## Relevant bibliography (as printed; the in-system, moon-tour and resonance references in full)
[3] A. Simon, F. Nimmo, and R. C. Anderson, "Planetary mission concept study for the 2023-2032 decadal survey:
Uranus Orbiter and Probe," National Aeronautics and Space Administration, Tech. Rep., 2021. [Online]. Available:
https://www.researchgate.net/publication/381251945
[4] L. N. Fletcher et al., "Ice giant system exploration within ESA's voyage 2050," Exp. Astron., vol. 54, pp.
1-11, 2022.
[5] R. A. Jacobson and R. S. Park, "The orbits of uranus, its satellites and rings, the gravity field of the
Uranian system, and the orientation of the poles of Uranus and its satellites," Astron. J., vol. 169, no. 2,
2025, Art. no. 65.
[7] H. Yang, J. Hu, X. Bai, and S. Li, "Review of trajectory design and optimization for Jovian system
exploration," Space: Sci. Technol., vol. 3, 2023, Art. no. 0036.
[11] S. J. Saikia et al., "Aerocapture assessment for NASA ice giants pre-decadal survey mission study," J.
Spacecraft Rockets, vol. 58, no. 2, pp. 505-515, 2021.
[12] A. P. Girija, "A flagship-class uranus orbiter and probe mission concept using aerocapture," Acta
Astronautica, vol. 202, pp. 104-118, 2023.
[13] D. C. Gochenaur, M. P. Jones, J. J. Norheim, and O. L. de Weck, "Orbit plane rotation using aerocapture,"
J. Spacecraft Rockets, pp. 1-11, 2025.
[15] X. Li, D. Qiao, and C. Circi, "Mars high orbit capture using manifolds in the Sun-Mars system," J.
Guidance, Control, Dyn., vol. 43, no. 7, pp. 1383-1392, 2020.
[16] F. Topputo and E. Belbruno, "Earth-Mars transfers with ballistic capture," Celestial Mechan. Dynamical
Astron., vol. 121, no. 4, pp. 329-346, 2015.
[23] N. Zhang, H. Baoyin, and H. Li, "Baseline trajectory design for solar polar detection using gravity
assists," Acta Astronautica, vol. 201, pp. 353-363, 2022.
[24] N. Zhang, D. Wu, and H. Baoyin, "Multiple ballistic resonant gravity-assist trajectory optimization using A*
algorithm with pruning," J. Guid., Control, Dyn., under review, 2025.
[25] A. Huang, H. Li, and Y. Luo, "Resonance selection and path planning for noncoplanar transfers via multiple
gravity assists," Adv. Space Res., vol. 73, no. 8, pp. 4213-4225, 2024.
[26] L. A. D'Amario, L. E. Bright, and A. A. Wolf, "Galileo trajectory design," Space Sci. Rev., vol. 60, pp.
23-78, 1992.
[27] D. L. Matson, L. J. Spilker, and J.-P. Lebreton, "The Cassini/Huygens mission to the Saturnian system,"
Space Sci. Rev., vol. 104, no. 1, pp. 1-58, 2002.
[28] M. Gavira-Aladro and C. Bombardelli, "Lambert-free solution of multiple-gravity-assist optimization
problem," J. Guid., Control, Dyn., vol. 47, no. 9, pp. 1822-1838, 2024.
[30] N. Zhang, W. Di, and H. Baoyin, "Analytic gradient computation and applications in multitarget multi-impulse
trajectory optimization," J. Guid., Control, Dyn., vol. 47, no. 2, pp. 347-357, 2024.
The remaining references ([1], [2], [6], [8]-[10], [14], [17]-[22], [29], [31]-[33]) are general (Uranus science,
Juno, Tianwen-1, low-thrust shaping, PSO, NLopt, Battin) and are not in-system Uranus or Neptune trajectory
works. There is no Neptune reference and no Triton content. Note that ref. [24], the A* resonant-GA method used
here, is printed as "under review, 2025".

## Suggested topology label
**mga-tour** (weak fit). Justification from the text: the Titania strategy is "multiple resonant GAs" of a single
moon to rotate the excess-velocity vector and change inclination during capture, a one-moon non-repeating
gravity-assist sequence rather than a repeating two-moon cycler, a V-infinity-leveraging pump tour, a resonant
periodic-orbit connection or a libration-orbit transfer. "None" would also be defensible, since the scope is
orbit capture rather than a moon tour.

## Provenance of this digest
Full text read from `pdftotext -layout` and the rendered pages 6-9 (Tables I-IV, Figures 6-9). Not checked: the
published version beyond this PDF, the contents of the cited UOP report [3].
