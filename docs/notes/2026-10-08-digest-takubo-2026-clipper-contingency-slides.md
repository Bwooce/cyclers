# Digest: Takubo, Campagnola, Pellegrini & Anderson, slide deck "Preliminary Contingency Trajectory Planning for Europa Clipper's Galilean Moon Tour" (JPL CL#26-0101) (#960, #943)

The presentation of the paper digested in `digest-takubo-2026-scitech.md` (AIAA SciTech 2026-1262, doi
10.2514/6.2026-1262). The slide footers carry the date "1/11/2026" and the export-control line "This document
has been reviewed for export control and it does NOT contain controlled technical data". The deck has no DOI
of its own.
- Source file: upload `db6ff460-CL26_0101.pdf`, 37 slides, md5 b9726c202ef832dd9a14c80e5878f409. It was made
  with Aspose.Slides for .NET 23.5 from PowerPoint. The text layer covers titles and bullets only. Equations,
  tables and plot labels are raster images.
- Proposed filename:
  `cyclers_pdf/papers/takubo-campagnola-pellegrini-anderson-2026-preliminary-contingency-trajectory-planning-europa-clipper-galilean-moon-tour-aiaa-2026-1262-doi-10.2514-6.2026-1262-slides-jpl-cl-26-0101.pdf`.
- How I read it:
  - all 37 slides rendered at 80 dpi and read on the images;
  - slides 4, 31 and 36 also at 220 dpi, cropped and zoomed;
  - each slide compared against the paper.
- Collision verdict: in `collision.md`. The deck adds no G-C or G-E sequence. **No collision with gc-1,
  gc-2, the k = 4-6 class or ge-1..3.**

## 0. Verdict

Most of the deck repeats the paper's figures: Figs. 2, 3, 4, 5, 6, 7, 8 and 9 and Table 3. **What it adds is
small but real:**
1. a different Star input table (slide 36);
2. COSMIC tolerances and a lockfile name (slide 37);
3. one-line verdicts on two other JPL tools (slide 27);
4. "Key takeaways" on the V_inf maps (slides 24-26);
5. a few numbers that differ from the paper (section 2).

Nothing here is a gate input. Use it as a companion file next to the paper. No catalogue or code proposal.

## 1. Content found only in the slides

| Slide | Content (read on the image) |
|---|---|
| 2-3 | Moon periods: "Callisto (T=16.7 days), Ganymede (T=7.15 days), Europa (T=3.55days), Io (T=1.77 days)"; "Almost coplanar & circular". A radiation-belt rendering is shown. |
| 4 | Cruise sketch: "LAUNCH NET 10/10/2024", Mars and Earth gravity assists, Jupiter orbit insertion; "Interplanetary transfer: 5.5 years". The gravity-assist dates are too low-resolution in the source raster to read reliably at 220 dpi; not transcribed. |
| 5 | The 21F31 tour plot (ref. [1] = Campagnola et al. 2025 JAS) with phase labels (COT-1..5, Petal Rotation, Sub-Jovian coverage, Disposal transfer, TEO1/Pump Down), and a radiation sketch from Buffington 2014 (ref. [2], held). |
| 6 | "Contingency Playbook"; "Potential failures: Radiation-related, Missed thrust etc." |
| 8 | Key questions: minimum Delta-V to escape; the Delta-V vs TID trade ("Can we reduce TID by adding Delta-V?"); escape from inclined orbits, "crank-over-the-top (COT): up to 5 deg. w.r.t. Jovian ecliptic plane". The paper gives no inclination number. |
| 13 | Star described as "Patched conics of pre-generated Lambert-like arcs – polynomial increase in search", with Lambert arcs (two-burn) and VILT; "Too restrictive -> no solutions; Too relaxing -> computation time explodes". |
| 14, 36 | "Important constraints" / "Custom-made constraints": "Conjunction: Earth-Sun-Jupiter angle", "Flyby strictly not allowed" (read as: no flyby during conjunction, INFERRED from the indentation); "Delta-V_lev < 5 m/s (and hope full-ephemeris correction can solve this)"; "Fix the inclination and epoch (dep./arr.) of the first leg"; "Bound RAAN of the first leg. When inclined, two possible RAANs (ascending/descending)". The paper has none of these. |
| 16 | COSMIC: the radiation is a "High-fidelity model based on the ephemeris". |
| 21 | Conclusion: "(Semi-)Ballistic escape is feasible if the inclination is zero; if not, it will likely spend < 100m/s of Delta-V." Future work: untargeted flybys, and "Consideration of eclipse at the patched-conics model" (eclipse is not in the paper). |
| 24 | Callisto V_inf map, "Key Takeaways": "4 flybys should be enough to reach the terminal orbit"; "2nd orbits with the highest rho: Eu 4:1 -> Ca 6:5; Eu 5:1 -> Ca 3:2; Eu 6:1 -> (Ca 3:2) -> Ca 5:3". |
| 25 | Ganymede V_inf map, "Key Takeaways": "Ganymede may contribute to the escape, especially if the V_inf is relatively low, but generally the energy w.r.t. Ganymede is too high to effectively change the V_inf w.r.t. other moons." |
| 26, 35 | Callisto V_inf map with annotations: Eu 4:1, Eu 5:1, Eu 6:1, Ga 3:1 and Ga 3:1+ markers; V_inf,Eu = 4 km/s and V_inf,Ga = 7 km/s contours; TID contours; the reachable-set curves (the "V_inf-contour of other moons" overlay). |
| 27 | "Previous methodologies": "STORM (Clipper's tour design tool): cannot accommodate large leveraging maneuvers"; "Dyno (Saturnian-moon tour optimization tool): would struggle with inter-moon transfer (focus on the ephemeris-free same-moon transfer)". |
| 28-31 | The Tisserand-graph build-up: TID per orbit (0-180 kRad colour scale), then the admissible domains, then the V_inf contours with maximum bending, then the resonances. Slide 31 marks the three nominal orbits with stars and the red V_inf,Ca contours (section 2 item 1). "* Smaller SC rev = less TID". |
| 32-35 | The V_inf-rho map build-up. Slide 33 shows the reachable region (red fill) and the bending formula "delta = 2 sin^-1(GM_fb / (GM_fb + r_p V_inf^2))". The paper's eq. (1b) has the same form. Slide 33's footer says "Yuji Takubo (392M) Final Presentation", not the export line. |
| 37 | COSMIC setup: "Multiple-shooter + correction process via SNOPT"; "lockfile: jup387_V6"; "Dynamics: N-body model (full-ephemeris), point mass, no SRP"; "postol = 1e-1 km, veltol = 1e-4 km/s (preliminary design)"; "No additional constraints". Pipeline: "Star solution file (.mat) -> COSMIC input file (.py)"; extract V_inf,arr, V_inf,dep and encounter/burn times; "Add CPs (periapsis in the multi-rev arc) and BPs (Delta-V = 0 at apoapsis), if needed"; "* Zoso (Star-compatible optimizer) struggled finding a feasible solution… could have been a good buffer?" |

## 2. Differences from the paper (numbers and settings)

1. **V_inf,Ca of the nominal orbits** (slide 31, image):
   - The slide prints "Eu 4:1 @V_inf,Eu = 4 km/s ~= @V_inf,Ca = 5 km/s", "Eu 5:1 ... ~= 5.7 km/s" and
     "Eu 6:1 @4.5 km/s ~= 6.3 km/s".
   - Paper Table 2 (p.5): 5.077, 5.825 and 6.376. Our circular coplanar check reproduces the paper's values
     exactly.
   - 6.3 and 5 are truncations of the paper's values, but 5.7 is not a rounding of 5.825. Treat the slide
     numbers as approximate. Use the paper.
2. **Star input table** (slide 36 vs paper Table 3, p.11):
   - N_rev[1]: the slide has "±1, only 1 rev, long & short". The paper has "(0, ..., ±m), up to m rev (from
     the nominal tour design)".
   - N_rev[2-4]: the slide has "(0, ±1, ±2, ±3), Up to 3 revs, long & short". The paper has (0, ±1, ..., ±5)
     if rho_Eu <= 5.5 and (0, ±1, ..., ±3) otherwise.
   - N_rev[5]: on the slide the two rows (0..±5) and (0..±3) and the rho <= 5.5 condition sit beside
     N_rev[5]. In the paper, N_rev[5] is (0, ±1, ..., ±3) only, and the split belongs to N_rev[2-4]. The
     slide layout may be misaligned. Use the paper.
   - Null(3): the slide comment is "B[2] -> B[4] is allowed". The paper says "skipping B[3] is allowed".
     They mean the same thing.
   - Comments: the slide says "from COSMIC inputs" where the paper says "from the nominal tour design".
     H_leg: the slide has "+Z, Only direct Z = [0,0,1]^T"; the paper has "+[0,0,1], only prograde orbits".
     Delta-V_enc: the slide comment is "Low Delta-V tolerance for nearly ballistic".
   - Red boxes on the slide mark alt_min[1] = alt_0, the 60 km Ganymede and Callisto minimum altitudes, and
     Delta-V_total = 0.200. These are presumably the values changed for the COT case (INFERRED).
   - All bounds are the same in both: V_inf ranges, TOF ranges, the 9e5 km r_min, the 0.01 steps and the
     crank grid.
3. **Safe-orbit set:**
   - Slide 12 gives the "Area of interest" as "Ca 5:3 / 2:1 / 5:2 / 3:1 (/ 7:2 / 4:1)" and the basis as
     "Perijove radius (~Callisto)".
   - The paper (p.7) says "Ca 5:3 to Ca 3:1" with perijove "the vicinity of Ganymede's orbital radius".
     The "~Callisto" on the slide looks like a slip.
4. **Eu 6:1 escape path:** slide 24 says "(Ca 3:2) -> Ca 5:3". The paper (p.9) gives "Ca 5:3, Ca 2:1, and
   Ca 3:1".
5. **Table 4 values** (slides 17-19): the same 18 numbers as the paper, with trailing zeros dropped (55.3,
   111.8, 189.5). Slide 18 is titled "Case 2: E19 (Petal, 5:1) – low TID with Delta-V".
6. **Moon periods** (slide 2): 16.7 / 7.15 / 3.55 / 1.77 d. Our values are 16.691 / 7.155 / 3.552 d for
   Callisto, Ganymede and Europa, so the slide values are consistent.

## 3. Citation mining

The deck cites three works:
- [1] Campagnola et al. 2025 JAS, "Europa Clipper Mission Design: Design of the 21F31 Reference Tour". Not
  held. Its conference form (ISSFD 2024) is held. It is wanted-list row 11, for attribution only.
- [2] Buffington 2014 AIAA 2014-4105. Held.
- [3] Landau, Campagnola & Pellegrini 2022 JAS, Star. Not held, not listed (see the paper digest).

The tools STORM, Dyno and Zoso are named without a citation. Nothing new for the wanted list beyond the
paper digest's suggestions.

*Filed as `cyclers_pdf/papers/takubo-campagnola-pellegrini-anderson-2026-preliminary-contingency-trajectory-planning-europa-clipper-galilean-moon-tour-aiaa-2026-1262-doi-10.2514-6.2026-1262-slides-jpl-cl-26-0101.pdf`.*
