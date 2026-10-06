# Digest: Schwaniger 1963, "Trajectories in the Earth-Moon Space with Symmetrical Free Return Properties" (#960 batch 30)

A. J. Schwaniger (Future Projects Branch, Aeroballistics Division, NASA Marshall SFC), NASA TN D-1833,
Lunar Flight Study Series Vol. 5, June 1963, NTRS 19630007117. No DOI. 32 PDF pages: 7 pages of text,
15 figures and 1 page of references.
- Filed as `cyclers_pdf/papers/schwaniger-1963-trajectories-earth-moon-space-symmetrical-free-return-properties-nasa-tn-d-1833-ntrs-19630007117.pdf`.
  - OCR copy (filed): `ntrs_19630007117-ocr.pdf`, md5
    109d8005e502b31f986252b0baeb5b58.
  - Original NTRS scan (supplied file `ntrs_19630007117.pdf`): md5 e1037b53344bf11699ad09068247a63d, 32 pp.
- How I read it:
  - Text: the OCR layer, all of it. The text is clean apart from some scrambled line breaks on p.3.
  - The figure pages are rotated, and their OCR is useless.
  - I read these figures on 200 dpi page images: Fig. 5 (p.12), Fig. 6 (p.13), Fig. 9 (p.16),
    Fig. 11 (p.18) and Fig. 15 (p.22).
  - Figs. 7, 8, 10 and 12-14 were not image-read. Their numbers below come from the text. The reference list (p.23) was read on the image.
- **There are no tables.** Every number is a text statement or a graph reading.
- Wanted-list row 26 ("Earth-Moon cycler ancestors; `#948` R4"). Removed in batch 30.

- Scripts and outputs:
  - `schwaniger1963_freereturn.py` reproduces the coplanar free returns (positive control). Output:
    `schwaniger1963_freereturn_output.txt`. The run was stopped after the 1938 km cases.
  - `schwaniger1963_periodic.py` locates and verifies the periodic orbit. Output:
    `schwaniger1963_periodic_output.txt`.
  - Two first attempts found nothing because they were set up wrongly, not because the orbit is
    absent. They are kept for the record:
    - `schwaniger1963_check.py` shoots from perigee. Its speed window was too narrow and its events too
      strict, and the only roots it found are far-side circumlunar orbits with 54 d periods.
    - `schwaniger1963_moonshoot.py` shoots from the Moon and found no root.

## 0. Verdict

**This is an Apollo-era parametric study of planar and 3-D symmetric free returns in the circular
restricted problem. It is also the source of one Earth-Moon cycler: a cislunar, counter-rotation
(retrograde) symmetric periodic orbit. That orbit passes 177 km above the Earth and about 460 km
above the Moon in every period, about once a month. I reproduced it.**

- **Model.**
  - The circular restricted three-body problem, integrated by Cowell's method.
  - Injection and re-entry are at perigee, r = 6555 km (100 n.mi. altitude), horizontal.
  - **The report states no mass ratio, distance or GM.** I used mu = 0.012150, L = 384,400 km and
    GM_E = 398,600.4 km^3/s^2. I also ran mu = 1/81.30 to test the sensitivity.
- **Two kinds of free return** (p.2), both from Miele's (1960) image theorem:
  - First kind: a perpendicular crossing of the Earth-Moon line at periselenum. The return mirrors
    the outbound leg about the x axis.
  - Second kind: a perpendicular crossing of the xz plane at periselenum. The return mirrors about
    the xz plane.
  - Coplanar orbits belong to both kinds.
  - "Circumlunar" means periselenum on the far side of the Moon. "Cislunar" means periselenum on the
    Earth side, after the path has gone outside the lunar orbit (Figs. 5, 6).
- **The periodic orbit (Sec. III.E, p.6-7).**
  - The paper argues: a coplanar free return is periodic if its injection is also on the x axis, at
    longitude 0 or 180 deg.
  - Cross-plotting Fig. 15 gives "one cislunar case ... counter-rotational injection at 180 degrees
    longitude and periselenum radius of about 2150 km", with a period "about 650 hours or about one
    month" (read off Fig. 11).
  - Summary (p.7): "a periodic trajectory which approaches the moon at about 2150 km radius
    approximately once each month."
- **Reproduction (mu = 0.012150, L = 384,400 km): the claim holds.**
  - Method: follow the counter-rotation cislunar free return (perigee 6555 km) in periselenum
    radius, and solve for the radius at which the perigee lies on the Earth-Moon line.
  - Result:
    - periselenum 2202.5 km on the Earth side of the Moon (464 km altitude);
    - perigee 6555.0 km at MEP longitude 180 deg (177 km altitude);
    - period 625.5 h = 26.06 d = 6.0015 time units;
    - Jacobi constant C = 1.08542 (the Jacobi drift over one period is 4e-11). This C omits the constant
      mu(1-mu); with it, C = 1.09727. The filer re-integrated the 6-digit state independently and found the
      half-period perpendicular crossing at t = 3.0008, perigee 6,560 km;
    - farthest distance 1.382 L (531,000 km) from Earth.
  - At the perigee, xdot = 1.5e-12, so the crossing is perpendicular. One full period closes to below
    1 m and 1e-4 m/s.
  - Initial condition at periselenum: x = 0.982120, y = 0, xdot = 0, ydot = -2.471279 (time unit
    375,190 s, speed unit 1.02455 km/s). The script holds the full-precision values.
  - Perigee speed: 10,960.0 m/s inertial, 10,977.4 m/s in the rotating frame. The inertial angular
    momentum at perigee is negative, so the orbit is retrograde ("counter-rotation", as printed).
  - With mu = 1/81.30 the result barely changes: 2229.8 km and 625.4 h.
  - Against the paper: the periselenum radius is 2.4% larger and the period 3.8% shorter. Both
    printed values are "about" figures read off graphs. My free-return times are 647.7 h at 2100 km
    and 626.0 h at 2200 km, which bracket the paper's 650 h near its 2150 km.
  - The rotating-frame shape (two large lobes crossing behind Earth, plus a small loop round Earth)
    matches the counter-rotation path drawn to scale in Fig. 6.
  - Stability: strongly unstable. A finite-difference monodromy gives a largest eigenvalue of about
    5e2. The estimate is rough: the unit pair is not recovered, because of the 177 km perigee.
- **Positive control (same model), coplanar free returns of the first kind, periselenum 1938 km.**

  | case | my round trip | paper | my injection longitude | Fig. 15 (image) |
  |---|---|---|---|---|
  | circumlunar, co-rotation | 138.5 h | 138-140 h (p.7) | - | - |
  | circumlunar, counter-rotation | 137.0 h | 138-140 h | - | - |
  | cislunar, co-rotation | 638.5 h | 620-700 h (p.7) | 168.1 deg | about 168 deg |
  | cislunar, counter-rotation | 688.9 h | 620-700 h | 197.9 deg | about 197 deg |

  The injection longitude is 360 deg minus the computed return longitude, by the first-kind mirror
  (p.2).
- **Unresolved: Fig. 9 injection speeds are about 64 m/s below mine.**
  - Fig. 9 (image) gives, at 1938 km, about 10,916 m/s (counter-rotation) and about 10,880 m/s
    (co-rotation). My rotating-frame values are 10,981.6 and 10,943.3 m/s.
  - The 36-38 m/s split between counter- and co-rotation matches only if Fig. 9 gives rotating-frame
    (MEP) speeds, since 2 x omega x r = 35 m/s. A constant offset of about 64 m/s then remains.
  - Times and longitudes match, so the offset is a speed convention or constant that the report does
    not state. It is not a different orbit.
- **Earth-Moon cycler: YES.**
  - It is periodic in the rotating frame and symmetric about the Earth-Moon line.
  - Every period it has one Earth perigee at 177 km altitude and one lunar pass at about 460 km
    altitude, well inside the lunar Hill radius. Its perigee of 6,555 km is 3 km below the lower edge
    of Vaquero's LEO-GEO insertion band of 6,558-42,164 km (catalogue line 56634).
  - It meets the brief's definition. It also meets the catalogue's working test for an Earth-Moon cycler, which needs
    the periselene inside the lunar SOI of 66,183 km (catalogue line 55649, Casoliva 1-2c, where the
    test fails).
- **Catalogue.**
  - No hit for "schwaniger", "D-1833" or "19630007117" in `data/catalogue.yaml`.
  - `arenstorf-em-figure8-1963` (line 9077) has null `jacobi_constant`, `period_nd` and `state_nd`
    (data_gaps from line 9154). This orbit should NOT be used to fill them. It is a retrograde
    cislunar free return with T = 6.00 TU and C = 1.085. The orbit usually shown as "the Arenstorf
    orbit" (Cook 2020, cited in the row) is a different orbit with T = 17.07 TU at mu = 0.012277.
  - No catalogue Earth-Moon row has C near 1.085 with T near 6.0.
    - The nearest Jacobi values are Casoliva 7-3b/7-3c (C = 1.0688, lines 56378 and 56483). Those
      are 7:3 resonant cyclers with seven perigees per period.
    - The period-6.0 rows (Braik-Ross line 55009, Vaquero lines 56730 and 56835) have C = 2.46-3.13.
    - Casoliva 2-1b (line 56063; C = 1.1964, perigee 6,790 km, periselene 92,590 km per the Liang 2020
      digest) is a prograde 2:1 resonant orbit with no lunar encounter. It is not a match.
    - Catalogue grep inventory: "earth-moon", "casoliva" and "arenstorf" hit only the Earth-Moon rows.
      "em-" gives 243 hits, mostly Earth-Mars rows (mcconaghy, aldrin-classic, s1l1, jones). The
      Earth-Moon rows are: 9077 Arenstorf; 9235 Genova-Aldrin 3-petal; 9469 Wittal 2022; 48219-48756
      Ross/Roberts-Tsoukkas; 48445 spatial (2,1); 48901-49164 Braik-Ross; 55643-56483 Casoliva;
      56588-57104 Vaquero. There are no hits for "schwaniger" or "hoelker".
    - Genova-Aldrin (line 9235) also has a 26-day lunar encounter period. It is a 3:1 prograde orbit
      with a 3,000 km perigee altitude, needs the Sun and maneuvers, and is labelled bicircular (line
      9239). Not a match.
    - Held digests checked for a collision ("retrograde", "free return", C near 1.0-1.1, periods near
      26 d or 6.0 TU, perigees near 6,5xx km):
      - Liang-Xu-Xu 2017: no hits.
      - Liang 2020: its nearest orbit has a 13,872 km perigee.
      - Casoliva 2008 seeds: 52a has C_J = 1.04619, but at mu = 1e-6, with no Moon passage data.
      - Genova-Aldrin mining note and the Arenstorf AIAA J digest: the 3:1 and Arenstorf m/k
        constructions are prograde.
      - None is this orbit.
    - So this is no duplicate.
  - **PROPOSAL only:** add one row,
    `schwaniger-1963-em-cislunar-retrograde-periodic-free-return`. Fields:
    - orbit_class cycler, model cr3bp, sense retrograde;
    - our_status known-reproduction (the orbit is published, 1963);
    - V0, not V1. The orbit itself closes (internal consistency). But the source prints only
      graph-read "about" values, and my reproduction matches them to 2.4% (radius) and 3.8% (period).
      That is outside the 1e-2 like-for-like tolerance the CR3BP rows use for V1 (line 55652). It is
      also strongly unstable, so a V2 multi-lap check would need a corrector, not plain propagation;
    - a `data_gaps` entry for the unstated mass ratio and the Fig. 9 speed offset.
  - It needs the usual literature-novelty and ratchet passes. It is NOT novel.
- **For `#948` R4:** this is a published, ballistic, Earth-grazing (177 km) and lunar-grazing
  (about 460 km) symmetric periodic orbit at the real mass ratio. It is a ready positive control for
  any both-primary-regularised corrector. The family continues in perigee radius: the paper's 6555 km
  is a design choice, not a constraint. Continued towards zero perigee radius, it should reach the
  double-collision orbits that R4 wants. That is my inference; the paper does not say it.

## 1. Content (READ)

- **3-D results (text, Figs. 7, 8, 10 and 12-14 not image-read).**
  - First kind: the inclination of the flight plane to the lunar MEP equator at periselenum is at most
    about 10.6-10.8 deg for any periselenum radius. It is 169.4 deg retrograde (circumlunar) or
    10.8 deg (cislunar) at 1938 km, and changes by less than 0.5 deg out to 20,000 km.
  - Second kind: the inclination grows with periselenum radius. It is 14 deg (cislunar) or 166 deg
    (circumlunar) at 1938 km, and any inclination from 0 to 180 deg is possible above 21,150 km.
- **Injection speed** (Fig. 9, image): about 10,850-10,917 m/s across 1,938-22,000 km, for both
  senses and both sides. It is nearly flat beyond 8,000 km. Read this with the offset caveat above.
- **Flight times** (Fig. 11, image; text p.7):
  - Circumlunar flights take 138-140 h near the Moon's surface. The time rises to about 220-235 h
    at 22 Mm.
  - Cislunar flights take 620-700 h near the surface. The time falls to about 260-265 h at 22 Mm.
  - The co- and counter-rotation curves cross near 3.4 Mm (cislunar) and near 3 Mm (circumlunar
    inset).
- **Injection loci** (Figs. 12-15): nearly circles for circumlunar flights and egg-shaped for
  cislunar ones. They move along the equator as the periselenum radius grows.
  - Fig. 15 (image), cislunar, second kind: the coplanar points are at about 168/197 deg (1,938 km),
    121/132 deg (3,738 km), about 88-91 deg (9,038 km) and about 72 deg (21,300 km).
  - The 1,938 km pair brackets 180 deg. This is the cross-plot the paper used for the periodic orbit.

## 2. Citation mining

- Miele, A. (1960), "Theorem of Image Trajectories in the Earth-Moon Space", Boeing Scientific
  Research Laboratories, Flight Sciences Laboratory Report No. 21, January 1960: not held, not on the
  wanted list.
  - It is the symmetry theorem behind both kinds of free return and behind every perpendicular-crossing
    corrector in the project.
  - New candidate, low priority: the result is textbook (Szebehely 1967, HELD). A journal version may
    exist; I have not verified one.
- Hoelker, Brand [printed so; probably Braud], Jean & Schwaniger (1960), "Parameter Sensitivity and Guidance Viewpoints for
  Circumlunar Flight and Return", MSFC MNN-M-AERO-3-60: not held. Not needed.
- Tucker 1962 (Lunar Flight Study Series Vol. 1, MTP-AERO-62-73), Schwaniger 1962 (Vol. 2,
  MTP-AERO-62-79), Braud 1962 (Vol. 3, MTP-AERO-62-84) and Braud 1963 (Vol. 4, MTP-AERO-63-7): none
  held, none wanted.
  - These are ephemeris-model mission studies. Not needed for the cycler lane.
- The paper does not cite Arenstorf 1963, Huang 1962 or Newton 1959. Hoelker & Winston 1968 (batch 30,
  same folder) does not cite this report either.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
