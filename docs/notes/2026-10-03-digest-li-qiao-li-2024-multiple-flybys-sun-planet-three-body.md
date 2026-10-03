# Digest — Li, Qiao & Li (2024), "Interplanetary Trajectory Design Using Multiple Flybys in the Sun-Planet Three-Body System" (Research Square preprint)

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; all 21 pages read; Table 3 header checked
against the rendered page because the text layer garbles it). **Citation:** Zhenyu Li, Dong Qiao, Xiangyu Li,
"Interplanetary Trajectory Design Using Multiple Flybys in the Sun-Planet Three-Body System",
Research Square preprint, posted 29 February 2024, **DOI 10.21203/rs.3.rs-3991174/v1**, School of
Aerospace Engineering, Beijing Institute of Technology (corresponding author D. Qiao). Licence CC BY 4.0.
Filed in the private paper corpus as
`li-qiao-li-2024-interplanetary-trajectory-design-multiple-flybys-sun-planet-three-body-system-research-square-preprint-doi-10.21203-rs.3.rs-3991174-v1.pdf`
(md5 `0c58728334a0d7b82dcf6831f9fdeaef`, 21 pages).

**Status caveat.** This is a Research Square preprint; the cover page carries no peer-review statement and the
paper itself does not say it has been peer reviewed. Treat it as not necessarily peer reviewed. No published
journal version is cited in the PDF.

**Project relevance in one line.** Adjacent, not a cycler paper. Every trajectory is a one-way Earth-to-outer-planet
transfer; nothing repeats (see question 4). Its value to us is as a mga-tour-class / Sun-Jupiter three-body
flyby-chaining method and as a source of a handful of printed transfer data.

## Q1. What is a "multiple flyby" and which types are found
- Dynamics model: a "patched planar restricted three-body problem (PCRTBP)" (Sec. 2.1). Each planet's system is a
  circular region xi of radius r_cr covering its satellites; the sphere-of-influence boundary is a circle of radius
  r_s = 3 r_soi, "where r_soi is the Hill radius of the planet". Inside a Sun-planet system S_i the motion is the
  planar CR3BP with normalised mass coefficient mu_i (Eq. 1) and Jacobi constant C_i (Eq. 2). Frames are converted
  between Sun-planet systems through DE440 (ref. [16]) via the heliocentric inertial frame.
- Definition (Sec. 2.2): "It is assumed that the spacecraft enters the planetary system from the outside of R and
  eventually escapes from this system." The spacecraft starts on the boundary P_e of region xi at t0 with state set
  by position angle alpha in [0, 2pi], velocity angle beta in [-pi, 0], Jacobi constant C_J and r_cr (Eqs. 3-4).
  The trajectory is propagated backward to t1 (reaching R) and forward to t2 (leaving R). "Before the spacecraft
  escapes the planetary system, there are p planet's pericenters within region xi", with constraints (Eq. 6)
  l >= r_s for t < t1, l = r_cr at t = t0, l >= r_s for t > t2, "and p should be larger than 1". So a multiple flyby is
  a single capture-like passage through the planet's neighbourhood containing p >= 2 pericentres before escape.
  (The text says p "larger than 1" for multiple flybys, yet also names a single-flyby type, so the p = 1 case is
  treated as a limiting case in the notation; not otherwise stated.)
- Notation: wp_{mn}^{p} (typeset as a script-p with p superscript, m, n subscripts): m = quadrant of entry, n = quadrant of
  escape, in the planet-centred rotating frame; p = number of flybys. Examples in Fig. 4, all Sun-Jupiter, C = 3.02:
  single-flyby wp_{32}^{1} (alpha = 108 deg, beta = -128 deg); two-flyby wp_{34}^{2} (alpha = 73 deg, beta = -121 deg);
  three-flyby wp_{32}^{3} (alpha = 52 deg, beta = -121 deg). (Degree signs are garbled in the text layer; values as
  extracted.) "The results in Fig. 3 show that two-flyby motions are easier to achieve than three-flyby motions."
- Two-flyby motions are split into eight families: wp_{32}^{2+}, wp_{32}^{2-}, wp_{34}^{2+}, wp_{34}^{2-},
  wp_{12}^{2+}, wp_{12}^{2-}, wp_{14}^{2+}, wp_{14}^{2-}. Superscript + = counter-clockwise flight relative to Jupiter at
  the first flyby, - = clockwise; a further superscript c = "the flight direction at the first flyby is the same as
  that at the second flyby", u = opposite direction (Fig. 5 a-h).
- Systems: the types are mapped for the **Sun-Jupiter** system only (normalised r_cr = 0.0041, r_soi = 0.0681,
  r_s = 0.2043; step sizes pi/180 for alpha and beta, 0.001 for C; normalised propagation times 2 forward and 5
  backward; at C = 3.02 there are 360 x 180 initial states, Fig. 3). Fig. 6 shows the fraction Theta of
  initial states giving two-flyby motion versus C: "the Jacobi constant C between 3.015 and 3.025 is more accessible to
  form two-flyby motions." The Neptune example also uses a **single-flyby** motion in the **Sun-Saturn** system
  (types not tabulated for Saturn).

## Q2. Trajectory patched conditions and connecting arcs between systems
- "alpha0, beta0, and C0, which can describe the initial state of multiple-flyby motions, are defined as trajectory
  patched conditions." They are found by brute-force gridding: any (alpha, beta, C) whose propagation satisfies Eq. 6
  is a patched condition (Sec. 3). These are pre-computed per system ("pre-designed by the method in Sec.3").
- Transfer assembly (Sec. 4.1): the usual patched-conic MGA timing constraints (Eq. 7, circular coplanar planets,
  heliocentric radius and phase match at each leg) are modified by adding a new design parameter, the in-system dwell
  time Dt_MF^i ("the flight time from the first pericenter to the last one of multiple-flyby trajectories", Eq. 8),
  which "relaxes" the phase constraint. Worked illustration: with a 3-year Dt_MF at Jupiter the spacecraft leaves
  Jupiter in 2037 instead of 2034, giving an intersection with Saturn's phase that a plain Jupiter gravity assist
  lacks (Fig. 8).
- Connection (Sec. 4.2): arcs in neighbouring systems are joined **by impulsive manoeuvres**, not ballistically. At the
  pericentre state x_per^i of system i, an impulse Dv_i^2 puts the spacecraft on the heliocentric-ish arc that reaches
  the patched entry state x_in^{i+1} of system i+1 (converted by coordinate transformation); a second impulse
  Dv_{i+1}^1 matches the velocity. Problem O1 (Eq. 9): minimise J = |Dv_i^2| + |Dv_{i+1}^1| subject to x(t_f) = x_in^{i+1}
  and Dt_{i+1} = JD_f - JD_d. The Earth-to-first-planet leg is designed in reverse, with Dv_1^1 at the first pericentre of
  the two-flyby trajectory so that the transfer is "tangent to the Earth's orbit in the sun-planet system" (Fig. 10 flow).
  The solver used for O1 is not stated.
- Not stated: the numerical integrator, tolerances, any Lambert or shooting detail, the Jupiter-Saturn heliocentric
  arc model between systems beyond "propagated forward tangent to Saturn's orbit in the Sun-Jupiter system".

## Q3. Transfers designed (Sec. 5)
Selection: wp_{34}^{2+} and wp_{34}^{2-} chosen "since the spacecraft's energy is expected to increase as it leaves the
Jupiter system [20]" (ref. [20] is Qi & Xu 2016 on lunar gravity assist in the restricted four-body problem);
Sun-Jupiter C_J in [3, 3.027].

**Earth-to-Jupiter (Table 1, transcribed exactly)**

| Types of transfer | alpha | beta | flight time (days) | Dv1 (km/s) | C |
|---|---|---|---|---|---|
| long-term | 163 deg | -67 deg | 6941.1409 | 1.4313 | 3.02 |
| short-term | 36 deg | -135 deg | 1102.4798 | 1.3960 | 3.02 |

wp_{34}^{2+} gives the long transfer, wp_{34}^{2-} the short one; "wp_{34}^{2-} is more suitable for designing the
Jupiter-Saturn trajectory, which is due to the shorter flight time." The paper states a conventional Jupiter-gravity-assist
Earth-Saturn route has no opportunities between 2030 and 2036 (Fig. 13, straight-line-escape argument); the method finds
four launch windows: "February 2030, March 2031, March 2032, and April 2033", with Dv < 1 km/s. Sequence: Earth launch,
Jupiter two-flyby (wp_{34}^{2-c} or wp_{34}^{2-u}), Saturn rendezvous. Two Dv: Dv1 at first Jupiter pericentre, Dv2 at last.

**Jupiter-Saturn opportunities (Table 2, transcribed exactly; dates as printed, month/day/year)**

| No. | Launch Date | Dv1 Date | Dv2 Date | Arrival Date |
|---|---|---|---|---|
| 1 | 2/05/2030 | 5/12/2032 | 3/20/2038 | 6/21/2047 |
| 2 | 2/04/2030 | 5/08/2032 | 3/13/2038 | 6/19/2047 |
| 3 | 2/28/2031 | 3/12/2033 | 4/06/2038 | 7/02/2047 |
| 4 | 3/02/2031 | 6/17/2033 | 3/09/2038 | 6/28/2047 |
| 5 | 3/23/2032 | 6/11/2034 | 2/07/2038 | 6/29/2047 |
| 6 | 3/24/2032 | 6/06/2034 | 1/25/2038 | 6/21/2047 |
| 7 | 4/21/2033 | 8/24/2035 | 12/21/2037 | 6/27/2047 |
| 8 | 4/23/2033 | 9/25/2035 | 12/12/2037 | 6/21/2047 |
| 9 | 3/01/2031 | 5/06/2033 | 10/21/2038 | 5/28/2047 |
| 10 | 3/02/2031 | 6/08/2033 | 7/27/2038 | 6/25/2047 |
| 11 | 3/25/2032 | 7/13/2034 | 10/7/2038 | 6/06/2047 |
| 12 | 3/24/2032 | 4/21/2034 | 7/16/2038 | 7/03/2047 |
| 13 | 4/20/2033 | 8/07/2035 | 1/04/2038 | 7/27/2047 |
| 14 | 4/19/2033 | 7/29/2035 | 1/13/2038 | 8/07/2047 |

**Flight time and impulses (Table 3, transcribed exactly).** Column headers as typeset: "Trajectory number", the family
symbol (wp_{12}^{2-} in the header cell, though every row is wp_{34}^{2-c} or wp_{34}^{2-u}; a likely typesetting
slip in the header, not asserted as an error), "T days", "Dt_PF days" (the text and Eq. 8 call this Dt_MF), "Dv1 km/s",
"Dv2 km/s", "Sum Dv km/s".

| No. | Type | T (days) | Dt_PF (days) | Dv1 (km/s) | Dv2 (km/s) | Sum Dv (km/s) |
|---|---|---|---|---|---|---|
| 1 | 2-c | 6345.1098 | 2138.1367 | 0.8019 | 0.0691 | 0.8710 |
| 2 | 2-c | 6343.7979 | 2135.0756 | 0.7909 | 0.0763 | 0.8672 |
| 3 | 2-c | 5967.6657 | 1850.8760 | 0.4552 | 0.0374 | 0.4926 |
| 4 | 2-c | 5962.0812 | 1725.1927 | 0.7089 | 0.0552 | 0.7641 |
| 5 | 2-c | 5576.3622 | 1336.4458 | 0.4519 | 0.0579 | 0.5098 |
| 6 | 2-c | 5567.0367 | 1329.2132 | 0.4125 | 0.0869 | 0.4995 |
| 7 | 2-c | 5180.0152 | 849.3724 | 0.3238 | 0.0610 | 0.3849 |
| 8 | 2-c | 5172.7083 | 808.9305 | 0.3979 | 0.0820 | 0.4799 |
| 9 | 2-u | 5932.4703 | 1994.4914 | 0.6841 | 0.0971 | 0.7813 |
| 10 | 2-u | 5960.6985 | 1874.7497 | 0.7499 | 0.0456 | 0.7956 |
| 11 | 2-u | 5551.7294 | 1547.4462 | 0.6901 | 0.1687 | 0.8587 |
| 12 | 2-u | 5578.8426 | 1546.9012 | 0.4036 | 0.0866 | 0.4903 |
| 13 | 2-u | 5210.2102 | 880.86765 | 0.4139 | 0.0479 | 0.4618 |
| 14 | 2-u | 5222.8855 | 899.02361 | 0.4381 | 0.0827 | 0.5209 |

(Type column: all rows are wp_{34}^{2-c} for 1-8 and wp_{34}^{2-u} for 9-14, order as in Table 2.) Text summary: flight time
"5172.7083 - 6345.1098 days, and Dv is 0.3849 - 0.8709 km/s"; the table's own largest sum is 0.8710, so the text's
upper bound differs from the table in the last digit (not an error we can resolve; both as printed).

**Worked Earth-Saturn example (Fig. 14, = Table 3 row 7 / Table 2 row 7).** Launch 21 April 2033 with v_inf at Earth
8.2355 km/s; 855.3733 days to the first Jupiter pericentre; Dv1 = 0.3238 km/s to enter the multiple-flyby trajectory
(C_J = 3.02, alpha = 103.25 deg, beta = -99.5 deg); 849.3724 days in the Jupiter system; Dv2 = 0.0610 km/s; 3475.2690
days to Saturn, arriving 27 June 2047 with v_inf at Saturn 1.1441 km/s. (855.3733 + 849.3724 + 3475.2690 = 5180.0147
days, consistent with T = 5180.0152 in Table 3 to rounding.)

**Earth-Jupiter-Saturn-Neptune (Sec. 5.2, Fig. 15-16, Table 4).** Jupiter two-flyby (wp_{34}^{2-}) then a **single-flyby** motion
in the Sun-Saturn system. Launch 23 May 2034, v_inf at Earth 8.5146 km/s; 855.3733 days to first Jupiter pericentre;
Dv1 = 0.3238 km/s; 849.3724 days; Dv2 = 0.0802 km/s into a "Jupiter-Saturn's SOI transfer trajectory"; 4795.0294 days
of flight, reaching the Saturn SOI target position on 8 March 2052 with v_inf at Saturn 2.3708 km/s. Then Dv3 = 0.825
km/s at the Saturn SOI boundary (to satisfy the single-flyby patched conditions; without it "the spacecraft cannot
naturally encounter Neptune"); 388.4188 days from the SOI to the Saturn pericentre; Dv4 = 0.7963 km/s at pericentre;
rendezvous with Neptune 2 January 2072, v_inf at Neptune 5.3697 km/s. "The total impulse of the Earth-to-Neptune
transfer ... is 2.0253 km/s, and the flight time is 37.6384 years." (0.3238 + 0.0802 + 0.825 + 0.7963 = 2.0253,
consistent.) Departure from LEO 200 km (r_pE about 6578 km): Dv_L = 6.1329 km/s (Eq. 10); capture into the orbit used
by Hughes et al. (ref. [21]), periapsis radius 61555 km, eccentricity 0.9737: Dv_C = 1.0362 km/s (Eq. 11).

Table 4 (transcribed exactly): "Multiple flybys", path JSN, Dv_L 6.13, Dv_C 1.04, Dv 2.03, Total Dv 9.19 km/s, launch
5/23/2034; "Multiple gravity assists", path JSN [21], Dv_L 7.78, Dv_C 2.26, Dv "-", Total Dv 10.04 km/s, launch 6/20/2058.
(The paper labels the path JSN, i.e. Jupiter-Saturn-Neptune.) Claim: the most recent ballistic MGA opportunity "is around
2058[21], which is advanced to 2034".

Observation worth recording (our reading, not a paper statement): the Neptune example reuses the Saturn example's exact
Jupiter-stage numbers (855.3733 d, Dv1 = 0.3238 km/s, 849.3724 d) although its launch date (2034) and Earth v_inf differ from
those of the 2033 Saturn example (the paper does not explain this; Dv2 does differ, 0.0802 vs 0.0610). Treat these as
printed values, with unexplained provenance, if ever used.

## Q4. Anything periodic, repeating, resonant or cycler-like?
**None.** Text search of the full paper: "cycler", "cyclic", "periodic", "resonan", "Tisserand", "Poincar", "invariant",
"return", "revolution", "Lambert" all occur zero times (the only "Lyapunov" is in a reference title, Canales et al.).
Every trajectory is a one-way Earth-to-target transfer with an impulsive-manoeuvre patch at each system boundary.
The deciding passages:
- Sec. 2.2: "It is assumed that the spacecraft enters the planetary system from the outside of R and eventually escapes
  from this system."
- Sec. 4.2: "The trajectory from Earth to the target planet is divided into n segments ... The first segment is the
  trajectory from Earth to the first planet, and the nth segment is the transfer from the (n - 1)th planet to the target planet."
- Sec. 6: "Interplanetary transfer is achieved by connecting multiple-flyby trajectories in different sun-planet systems."
The two or three pericentre passages are a single transient capture-and-escape episode in the Jupiter (or Saturn) region,
not a repeating encounter pattern; the paper never closes a trajectory back on itself, never revisits Earth, and never
mentions a repeat period. Its "opportunities" are launch windows (one-off), not a recurring itinerary. One tangential
remark: the sentence "spacecraft have the capability to enter the Jupiter system more than once" (Sec. 5.2) refers to
the in-system multiple flybys, not to a recurrent route.

## Q5. Reproduction targets
Yes, a small set, all one-way transfers and all under the paper's own model (planar PCRTBP patched with DE440
transformations, impulsive patch, unstated optimiser), so a reproduction needs the paper's coordinate-patching
conventions, which are only sketched (refs [17, 18, 19]). Printed data usable as targets: Table 1 (two Earth-Jupiter
rows), Tables 2 and 3 (14 Jupiter-Saturn rows with dates and Dv), the Saturn worked example (Fig. 14 numbers) and the
Neptune worked example with Table 4. Not printed: state vectors, ephemeris epochs of the pericentres to better than the
dates above, the Saturn single-flyby patched conditions (alpha, beta, C), r_cr for Saturn, any integrator tolerances,
and the Fig. 3 phase-space map in numeric form. The only three-body invariants stated numerically: Sun-Jupiter r_cr =
0.0041, r_soi = 0.0681, r_s = 0.2043 (normalised), C = 3.02 as the working Jacobi constant, and the 3.015-3.025 band for
two-flyby motion. Because these are not cycler data, none is a V0/V1 candidate for any current catalogue row.

## Q6. Bibliography (21 references, in full as printed)
1. Agarwal, S., Pandey, R., Jeyaseelan, C.: Exploring the effect of various plasma parameters on whistler mode growth rates in the jovian magnetosphere. Astrophysics and Space Science 364, 1-9 (2019) doi 10.1007/s10509-019-3623-z
2. Penner, A.R.: A proposed experiment to test gravitational anti-screening and mond using sun-gas giant saddle points. Astrophysics and Space Science 365, 1-16 (2020) doi 10.1007/s10509-020-03870-x
3. Atreya, S.K., Hofstadter, M.D., Reh, K.R., In, J.H.: Icy giant planet exploration: Are entry probes essential? Acta Astronautica 162, 266-274 (2019) doi 10.1016/j.actaastro.2019.06.020
4. Liu, X., Schmidt, J.: Dust in the jupiter system outside the rings. Astrodynamics 3, 17-29 (2019)
5. Sun, Y., Zhao, J., Hou, C., Jiao, W.: Highlight advances in planetary physics in the solar system: In situ detection over the past 20 years. Space: Science & Technology 3, 0007 (2023) doi 10.34133/space.0007
6. Yang, H., Hu, J., Bai, X., Li, S.: Review of trajectory design and optimization for jovian system exploration. Space: Science & Technology 3, 0036 (2023) doi 10.34133/space.0036
7. Kohlhase, C., Penzo, P.A.: Voyager mission description. Space Science Reviews 21(2), 77-101 (1977) doi 10.1007/BF00200846
8. Adriani, A., Moriconi, M., Mura, A., Tosi, F., Sindoni, G., Noschese, R., Cicchetti, A., Filacchione, G.: Juno's earth flyby: the jovian infrared auroral mapper preliminary results. Astrophysics and Space Science 361, 1-8 (2016) doi 10.1007/s10509-016-2842-9
9. Li, X., Qiao, D., Chen, H.: Interplanetary transfer optimization using cost function with variable coefficients. Astrodynamics 3, 173-188 (2019) doi 10.1007/s42064-018-0043-8
10. Vasile, M., De Pascale, P.: Preliminary design of multiple gravity-assist trajectories. Journal of Spacecraft and Rockets 43(4), 794-805 (2006) doi 10.2514/1.17413
11. Wagner, S., Wie, B.: Hybrid algorithm for multiple gravity-assist and impulsive delta-v maneuvers. Journal of Guidance, Control, and Dynamics 38(11), 2096-2107 (2015) doi 10.2514/1.G000874
12. Rudd, R., Hall, J., Spradlin, G.: The voyager interstellar mission. Acta Astronautica 40(2-8), 383-396 (1997) doi 10.1016/S0094-5765(97)00146-X
13. Englander, J.A., Conway, B.A., Williams, T.: Automated mission planning via evolutionary algorithms. Journal of Guidance, Control, and Dynamics 35(6), 1878-1887 (2012) doi 10.2514/1.54101
14. Gad, A., Abdelkhalik, O.: Hidden genes genetic algorithm for multi-gravity-assist trajectories optimization. Journal of Spacecraft and Rockets 48(4), 629-641 (2011) doi 10.2514/1.52642
15. Abdelkhalik, O., Gad, A.: Dynamic-size multiple populations genetic algorithm for multigravity-assist trajectory optimization. Journal of Guidance, Control, and Dynamics 35(2), 520-529 (2012) doi 10.2514/1.54330
16. Park, R.S., Folkner, W.M., Williams, J.G., Boggs, D.H.: The jpl planetary and lunar ephemerides de440 and de441. The Astronomical Journal 161(3), 105 (2021) doi 10.3847/1538-3881/abd414
17. Canales, D., Howell, K.C., Fantino, E.: Transfer design between neighborhoods of planetary moons in the circular restricted three-body problem: the moon-to-moon analytical transfer method. Celestial Mechanics and Dynamical Astronomy 133(8), 36 (2021) doi 10.1007/s10569-021-10031-x
18. Canales, D., Howell, K.C., Fantino, E., Gilliam, A.J.: Transfers between moons with escape and capture patterns via lyapunov exponent maps. Journal of Guidance, Control, and Dynamics 46(11), 2133-2149 (2023) doi 10.2514/1.G007195
19. Qi, Y., Ruiter, A.: Short-term capture of the earth-moon system. Monthly Notices of the Royal Astronomical Society 476(4), 5464-5478 (2018) doi 10.1093/mnras/sty665
20. Qi, Y., Xu, S.: Study of lunar gravity assist orbits in the restricted four-body problem. Celestial Mechanics and Dynamical Astronomy 125, 333-361 (2016) doi 10.1007/s10569-016-9686-z
21. Hughes, K.M., Moore, J.W., Longuski, J.M.: Preliminary analysis of ballistic trajectories to neptune via gravity assists from venus earth mars jupiter saturn and uranus. In: AAS/AIAA Astrodynamics Specialist Conference (2013)

Bibliography observations: the paper cites **no** Tisserand-Poincare, resonant-flyby, Lambert-cycler or Strange-Longuski
works (the only resonant-style machinery is implicit in [17-19]). The three-body flyby lineage it actually builds on is
[19] Qi & Ruiter 2018 (alpha/beta/C initial-state parametrisation of short-term capture; the definitions of alpha and
beta are delegated to it) and [20] Qi & Xu 2016 (energy change on leaving the secondary's system); [17, 18] supply the
rotating-inertial frame transformations; [21] Hughes-Moore-Longuski 2013 is the ballistic JSN MGA baseline for Table 4.
For our corpus-index check: Hughes-Moore-Longuski 2013 and Qi's 2016/2018 papers are the candidates to grep CORPUS_INDEX
for; the paper itself gives no corpus-relevant moon-tour or cycler citations. Reference [6] (Yang et al. 2023 Jovian
review) is already digested in this project.

## Status
Filed; digest written. Classification suggestion for the index owner (not applied here): adjacent literature,
mga-tour-class interplanetary method (Sun-Jupiter planar CR3BP multiple-pericentre flyby plus impulsive patching);
not a cycler, no catalogue row implied.
