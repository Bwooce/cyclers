# Digest: Strange & Longuski 2002, "Graphical Method for Gravity-Assist Trajectory Design" (#960 batch 26)

N. J. Strange and J. M. Longuski (Purdue), Journal of Spacecraft and Rockets 39(1):9-16 (2002), doi
10.2514/2.3800. Received 23 January 2001.
- Filed as `cyclers_pdf/papers/strange-longuski-2002-graphical-method-gravity-assist-trajectory-design-jsr-39-1-9-doi-10.2514-2.3800.pdf`.
  8 pages, text layer, md5 cfb8f8e43179be1f0b63dc4ec935c89b. Supplied by the owner.
- I read the method, the conclusions and the references from the text layer. Tables 1-3 (flight times)
  were not transcribed.
- Wanted list: it was in the Golubev-references row (Tier D). That row now drops it.

## 0. Verdict

**The origin paper for "Tisserand graphs" (the P-rp plot), as background for the gauntlet and the tour
bookkeeping.** No cycler content.
- A flyby rotates v_inf at fixed magnitude. So, relative to one planet, all reachable heliocentric orbits
  lie on one v_inf contour, which is a contour of constant Tisserand parameter.
  - The contours are plotted as period P against periapsis r_p, for circular coplanar planets.
  - The far right end of a contour has alpha = 0 (v_inf aligned with the planet's velocity: highest
    energy, flyby at perihelion). The left end has alpha = 180 deg (lowest energy, flyby at aphelion).
- **Tick marks** mark the maximum orbit change per flyby at a minimum altitude: 200 km at the
  terrestrial planets, 5 R_J at Jupiter. Counting tick spacings gives the minimum number of flybys of
  that body. For example, a VEEGA needs two Earth flybys to rotate a 9 km/s v_inf.
- **Intersections of contours of different planets are the possible transfer orbits between them.**
  Paths that connect on the graph are energy-feasible.
  - **"Paths that do not exist on a Tisserand graph are strictly infeasible for ballistic
    trajectories."** It is a necessary, not sufficient, condition, because phasing is ignored.
  - Phasing is left to path solvers (STOUR).
- Flight-time estimates (Tables 1-3) bound the shortest ballistic path to each planet for launch v_inf
  = 3, 5, 7, 9 and 11 km/s, ignoring phasing.
  - Resonant returns up to 6:1 are included for the inner planets, plus 3:2, 1:2 and 1:3 for Venus,
    Earth and Mars.
  - Paths are truncated after five bodies. Without a limit, 15- and 17-body paths to Pluto appear, with
    infeasible phasing.
- **Use in this project:**
  - The ballistic-feasibility screen: a cycler leg pair must lie on intersecting v_inf contours, with
    the turn bounded by the altitude tick spacing. This is the same physics as the project's Tisserand
    and turn-angle checks.
  - Cite this paper (with Labunsky et al. 1998 and Hollenbeck 1975 as earlier graphical forms) for that
    criterion.
- **Minovitch citation:** ref. 5 gives "JPL TR 32-464, Oct. 1963". So three sources say 32-464 and only
  Niehoff 1965 says 32-468.

## 1. Citation mining (refs 1-20)

- Held: Szebehely 1967 [19].
- Already on the wanted list (methods and textbook rows): Battin 1959 [3], Minovitch TR 32-464 [5] and
  Labunsky 1998 [18].
- Added in one Tier D row: Petropoulos, Longuski & Bonfiglio 2000 (JSR 37(6):776, Venus-Earth-Mars
  assists to Jupiter) [10], Longuski & Williams 1991 (CMDA 52:207, the automated gravity-assist design
  behind STOUR) [13], and Deerwester 1966 (JSR 3(10):1564, Jupiter swingby to the outer planets) [6].
  The wanted Deerwester row is his 1965 Venus-swingby paper, a different item.
- Not added (background or mission reports): Broucke AIAA 88-4220 [1], Roy [2], Sedov 1960 [4],
  Flandro 1966 [7], Farquhar & Stern 1990 [8], Hollenbeck AAS 75-087 [9], Rinderle JPL D-263 [11], the
  Purdue MS theses [12, 14, 15], Heaton et al. 2002 [16], Weinstein & Stetson AAS 89-433 [17], and
  Johnson & Longuski 2002 [20].
